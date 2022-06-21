#Haley Greenyer 
#Updated: 6/16/2022


#' Run Quality Control, make figures
#'
#' @param salmon_counts count matrix
#' @param dds DESeqDataSet object
#' @param vst VST count matrix 
#' @param annotation sample annotation df
#' @param grp_by attribute to group by, e.g. 'condition', 'replicate', etc.
#' @param out_dir option to specify directory to which figures and tables will be exported
#' @param cell option to specify cell line 
#'
#' @return QC figures and tables 
#' @export
#'
#' @examples run_QC(salmon_counts, dds, vst, annotation)
run_QC <- function(salmon_counts, dds, vst, annotation, grp_by = 'condition', out_dir=NULL, cell=NULL){
  
  vst_mat <- assay(vst)
  
  #barplot of library sizes
  librarySizes <- colSums(salmon_counts) 
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell, 'LibrarySizes.pdf',sep=''))
    barplot(librarySizes, 
            names=names(librarySizes), 
            las=2, 
            main="Barplot of library sizes")
    
    dev.off()
  }
  else{
    
    barplot(librarySizes, 
            names=names(librarySizes), 
            las=2, 
            main="Barplot of library sizes")
    
  }
  
  
  
  logcounts <- log2(vst_mat + 1)
  # make a color vector
  statusCol <- as.numeric(factor(annotation$condition)) + 1
  # Check distributions of samples using boxplots
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell,'SampleDistributions.pdf', sep=''))
    boxplot(logcounts, 
            xlab="", 
            ylab="Log2(counts)",
            las=2,
            col=statusCol)
    dev.off()
  }
  else{
    boxplot(logcounts, 
            xlab="", 
            ylab="Log2(counts)",
            las=2,
            col=statusCol)
  }
  
  
  GeneCounts <- counts(dds)
  idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
  sum(idx.nz)
  
  #estimate the size factors
  DESeq2Table <- estimateSizeFactors(dds)
  sizeFactors(DESeq2Table)
  #densities of counts
  multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],
                xlab="mean counts", xlim=c(0, 1000))
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell,'CountDensities.pdf', sep=''))
    multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],
                 xlab="mean counts", xlim=c(0, 1000))
    dev.off()
  }
  else{
    multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],
                 xlab="mean counts", xlim=c(0, 1000))
  }
  
  
  #empircal cumulative distribution functions (ECDFs), essential integrals of the densities and give the probability of
  #observing a certain number of counts equal to x or less given the data.
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell,'ECDFplot.pdf',sep=''))
    multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
               ylab="estimated probability(Fn(x))", xlab="mean counts", xlim=c(0, 1000))
    dev.off()
  }
  else{
    multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
               ylab="estimated probability(Fn(x))", xlab="mean counts", xlim=c(0, 1000))
  }
  
  
  
  #sample correlations
  d <- cor((vst_mat), method="spearman")
  hc <- hclust(dist(1-d))
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell,'SampleCorrelation.pdf',sep=''))
    ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)
    dev.off()
  }
  else{
    ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)
  }
  
  
  vst_cor <- cor(vst_mat)
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell,'CorHeatmap.pdf',sep=''))
    pheatmap(vst_cor)
    dev.off()
  }
  else{
    pheatmap(vst_cor)
  }
  
  
  #PCA
  pcaData <- plotPCA(vst, intgroup=c("condition", "replicate"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell,'PCA.pdf',sep=''))
    ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()
    dev.off()
  }
  else{
    ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()
  }
  
}

salmon_QC <- function(salmon_data, dds, vst, annotation, grp_by = 'condition', out_dir=NULL, cell=NULL) UseMethod('salmon_QC', salmon_data)

#' Salmon Quality Control from txi list 
#'
#' @param salmon_data txi list of salmon data 
#' @param dds DESeqDataSet object
#' @param vst VST DESeqDataSet
#' @param annotation sample annotation df
#' @param grp_by attribute to group by, e.g. 'condition', 'replicate', etc.
#' @param out_dir option to specify directory to which figures and tables will be exported
#' @param cell option to specify cell line 
#'
#' @return QC figures and tables
#' @export
#'
#' @examples salmon_QC(txi_list, dds, vst, annotation_df)
setMethod("salmon_QC",
          signature(salmon_data = "list"),
          function (salmon_data, dds, vst, annotation, grp_by = 'condition', out_dir=NULL, cell=NULL) 
          {
            run_QC(salmon_counts = salmon_data$counts, dds=dds, vst=vst, annotation=annotation, 
                   grp_by=grp_by, out_dir=out_dir, cell=cell)
          }
)


#' Salmon Quality Control from matrix
#'
#' @param salmon_data matrix of counts
#' @param dds DESeqDataSet object
#' @param vst VST DESeqDataSet
#' @param annotation sample annotation df
#' @param grp_by attribute to group by, e.g. 'condition', 'replicate', etc.
#' @param out_dir option to specify directory to which figures and tables will be exported
#' @param cell option to specify cell line 
#'
#' @return QC figures and tables
#' @export
#'
#' @examples salmon_QC(count_matrix, dds, vst, annotation_df)
setMethod("salmon_QC",
          signature(salmon_data = "matrix"),
          function (salmon_data, dds, vst, annotation, grp_by = 'condition', out_dir=NULL, cell=NULL) 
          {
            run_QC(salmon_counts = salmon_data, dds=dds, vst=vst, annotation=annotation, 
                   grp_by=grp_by, out_dir=out_dir, cell=cell)
          }
)

#' Salmon Quality Control from double  
#'
#' @param salmon_data matrix of counts
#' @param dds DESeqDataSet object
#' @param vst VST DESeqDataSet
#' @param annotation sample annotation df
#' @param grp_by attribute to group by, e.g. 'condition', 'replicate', etc.
#' @param out_dir option to specify directory to which figures and tables will be exported
#' @param cell option to specify cell line 
#'
#' @return QC figures and tables
#' @export
#'
#' @examples salmon_QC(counts, dds, vst, annotation_df)
setMethod("salmon_QC",
          signature(salmon_data = "double"),
          function (salmon_data, dds, vst, annotation, grp_by = 'condition', out_dir=NULL, cell=NULL) 
          {
            run_QC(salmon_counts = salmon_data, dds=dds, vst=vst, annotation=annotation, 
                   grp_by=grp_by, out_dir=out_dir, cell=cell)
          }
)

