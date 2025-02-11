% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/salmonQC.R
\name{run_QC}
\alias{run_QC}
\title{Run Quality Control, make figures}
\usage{
run_QC(
  salmon_counts,
  dds,
  vst,
  annotation,
  grp_by = "condition",
  out_dir = NULL,
  cell = NULL
)
}
\arguments{
\item{salmon_counts}{count matrix}

\item{dds}{DESeqDataSet object}

\item{vst}{VST count matrix}

\item{annotation}{sample annotation df}

\item{grp_by}{attribute to group by, e.g. 'condition', 'replicate', etc.}

\item{out_dir}{option to specify directory to which figures and tables will be exported}

\item{cell}{option to specify cell line}
}
\value{
QC figures and tables
}
\description{
Run Quality Control, make figures
}
\examples{
run_QC(salmon_counts, dds, vst, annotation)
}
