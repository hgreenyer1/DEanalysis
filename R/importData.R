#Haley Greenyer 
#Updated: 6/16/2022


#' Get Referemce
#'
#' @param gtf_dir Directory where GTF file is located  
#'
#' @return A txdb object that stores transcript metadata 
#' @export
#'
#' @examples get_ref('/path/to/gtf_file')
get_ref <- function(gtf_dir){
  
  txdb <- makeTxDbFromGFF(gtf_dir)
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
  return(tx2gene)
}


#' Make Annotation, import_salmon helper function
#'
#' @param salmon_dir Directory where salmon output file folders are located
#'
#' @return A data frame of sample annotations
#' @export
#'
#' @examples make_annotation('/path/to/salmon_folders')
make_annotation <- function(salmon_dir){
  
  sample <- list.files(salmon_dir) %>% str_replace_all( ".salmon_quant", "")
  cell = sapply(strsplit(sample, '_'), `[`, 1)
  condition = sapply(strsplit(sample, '_'), `[`, 2)
  replicate = sapply(strsplit(sample, '_'), `[`, 3)
  directory = list.files(salmon_dir)
  #make df
  annotation <- data.frame(sample = sample, directory= directory, cell = cell, condition = condition, replicate = replicate, stringsAsFactors = TRUE)
  rownames(annotation) <- annotation$sample
  return(annotation)
}


#' txi_salmon, import_salmon helper function
#'
#' @param salmon_dir Directory where salmon output file folders are located
#' @param annotation Data frame containing sample annotations
#' @param tx2gene A txdb object 
#'
#' @return a (txi) list containing matrices extracted from salmon output: abundance, counts, length
#' @export
#'
#' @examples txi_salmon('~/salmonDirectory',annotation_dir,tx2gene)
txi_salmon <- function(salmon_dir, annotation, tx2gene){
  files <- file.path(salmon_dir,annotation$directory, "quant.sf")
  names(files) <- paste0(annotation$sample)
  #all(file.exists(files))
  txi_salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
  return(txi_salmon)
}


#' Import Salmon Data
#'
#' @param salmon_dir Directory where salmon output file folders are located
#' @param tx2gene A txdb object 
#' @param out_dir option to specify output directory 
#'
#' @return a list containing two objects: an annotations df and corresponding salmon data extracted via txi_salmon
#' @export
#'
#' @examples import_salmon_data('~/salmonDirectory', tx2gene)
import_salmon_data <- function(salmon_dir, tx2gene, out_dir=NULL){
  annotation_df <- make_annotation(salmon_dir=salmon_dir)
  txi_salmon <- txi_salmon(salmon_dir = salmon_dir, annotation = annotation_df, tx2gene = tx2gene)
  
  #option to export raw counts 
  if(!is.null(out_dir)){
    file_path = paste(out_dir,'/','all_counts_raw.csv',sep='')
    write(txi_salmon$counts, file=file_path)
  }
  
  return(list('annotation'=annotation_df, 'txi_salmon'=txi_salmon))
}


   
#' Get Cell Subset
#'
#' @param cell name of cell line to subset
#' @param txi_salmon txi list
#' @param out_dir option to specify output directory 
#'
#' @return a list containing two objects: an annotations df and corresponding matrix of abundances
#' @export
#'
#' @examples get_cell_subset('Kasumi', txi_list)
get_cell_subset <- function(cell, txi_salmon, out_dir=NULL){
  #subsetting
  df.salmon_A = as.data.frame(txi_salmon$counts) 
  df.A = df.salmon_A %>% dplyr::select(contains(cell))
  df.A <- cbind(rownames(df.A), data.frame(df.A, row.names=NULL))
  names(df.A)[1] <- "GeneID"
  salmon.df <- df.A %>% 
    mutate(Gene = str_replace(GeneID, ".[0-9]+$", ""))
  salmon.df2 <- salmon.df[,c(26,2:25)] 
  mat.A <- as_tibble(salmon.df2) %>% as_matrix()
  
  #annotation
  file = colnames(mat.A)
  condition = sapply(strsplit(colnames(mat.A), '_'), `[`, 2) 
  replicate = sapply(strsplit(colnames(mat.A), '_'), `[`, 3) 
  #make df
  annotation <- data.frame(file = file, condition = condition, replicate = replicate, stringsAsFactors = TRUE)
  row.names(annotation$file)
  
  #option to export raw abundances -does not export annotations, that can be exported after QC
  if(!is.null(out_dir)){
    file_path = paste(out_dir,'/',cell,'Counts_raw.csv', sep='')
    write(mat.A, file=file_path)
  }
  

  return(list('annotation'=annotation, 'abundances'=mat.A))
  
}
